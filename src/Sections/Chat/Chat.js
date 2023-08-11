import React, { Component } from 'react';
import { useState, useEffect, createRef } from 'react';

import { Navigate, Link } from 'react-router-dom';

import Select from 'react-select';

import { variables, AG_GRID_LOCALE_RU } from '../Variables.js';

import update from 'immutability-helper';


export class Chat extends Component {

  constructor(props) {
    super(props);

    this.state = {
      token: variables.token,
      allow_page: variables.allow,
      loading: false,
      messages: [{'request': 'test 1', 'response': 'response 1', status: 200}, {'request': 'test 2', 'response': 'response 2', status: 200}],
      message: null,
      messageStatus: 0,
      sendMessage: null,
    }
  }

  getResponse = (task_id, interval = 1000) => {
    fetch(variables.API_URL + `/api/chat?task_id=${task_id}`,
      {
        headers: {
          'Content-Type': 'application/json;charset=utf-8',
          'Authorization': `Token ${variables.token}`,
        },
      }
    )
      .then((res) => {
        if (res.status == 202) {
          setTimeout(() => {
            return this.getResponse(task_id, interval)
          }, interval);
        } else if (res.status == 200) {
          return res.json()
        } else {
          throw Error(res.statusText)
        }
      })
      .then((data) => {
        try {
          this.setState({
            message: data.data,
            messageStatus: 200,
            loading: false,
          });
          this.state.messages.push({'request': this.state.sendMessage, 'response': data.data, status: 200})
        } catch {
          console.log('access')
        }
      })
      .catch((err) => {
        console.log(err);
        this.setState({ message: 'Произошла ошибка', loading: false, messageStatus: 500 });
        this.state.messages.push({'request': this.state.sendMessage, 'response': 'Произошла ошибка', status: 500})
      });
  }

  getRequest() {
    this.setState({ loading: true })
    fetch(variables.API_URL + '/api/chat', {
      method: 'POST',
      headers: {
        'Accept': 'application/json',
        'Content-Type': 'application/json;charset=utf-8',
        'Authorization': `Token ${variables.token}`,
      },
      body: JSON.stringify({
        message: this.state.sendMessage
      })
    })
      .then((res) => {
        console.log(res.status)
        if (res.ok) {
          return res.json()
        } else {
          throw Error(res.statusText)
        }
      })
      .then((result) => {
        var task_id = result.data;
        this.setState({ message: 'Ждем ответа', messageStatus: 201, loading: true })
        this.getResponse(task_id);
      })
      .catch((err) => {
        console.log(err);
        this.setState({
          message: 'Ошибка при отправке сообщения',
          messageStatus: 500,
          loading: false,
        });
      });
  }

  changeQueryText = (e) => {
    this.setState({ sendMessage: e.target.value });
  }

  updateMessageTape(message) {
    this.state.messages.push(message)
  }

  render() {
    const {
      token,
      loading,

      messages,
      messageStatus,
      message,
      sendMessage,

      allow_page,

    } = this.state;

    if (!token) {
      return <Navigate push to="/login" />
    } else if (allow_page === 0) {
      return <Navigate push to="/tematic_review" />
    } else if (allow_page === 1) {
      return <Navigate push to="/ddi_review" />
    } else {
      return (
        <>
          <header>
            <nav class="bg-white border-gray-200 px-4 lg:px-6 py-2.5">
              <div class="flex flex-wrap justify-between items-center">
                <div class="flex justify-start items-center">
                  <a href="" class="flex mr-4">
                    <img src="https://flowbite.s3.amazonaws.com/logo.svg" class="mr-3 h-8" alt="FlowBite Logo" />
                    <span class="self-center text-2xl font-semibold whitespace-nowrap">EBM Sechenov DataMed.AI</span>
                  </a>
                  {allow_page === 3 ?
                    <ul class="flex font-medium flex-row space-x-8">
                      <Link to="/tematic_review">
                        <li>
                          <a href="#" class="block py-2 pl-3 pr-4 text-gray-900 rounded hover:bg-gray-100 md:hover:bg-transparent md:hover:text-blue-700 md:p-0">Тематический анализ</a>
                        </li>
                      </Link>
                      <Link to="/chat">
                        <li>
                          <a href="#" class="block py-2 pl-3 pr-4 text-gray-900 bg-blue-700 rounded md:bg-transparent md:text-blue-700 md:p-0" aria-current="page">Поговорим</a>
                        </li>
                      </Link>
                      <Link to="/ddi_review">
                        <li>
                          <a href="#" class="block py-2 pl-3 pr-4 text-gray-900 rounded hover:bg-gray-100 md:hover:bg-transparent md:hover:text-blue-700 md:p-0">Факты для EBM</a>
                        </li>
                      </Link>
                      {variables.admin?
                        <Link to="/admin">
                          <li>
                            <a href="#" class="block py-2 pl-3 pr-4 text-gray-900 rounded hover:bg-gray-100 md:hover:bg-transparent md:hover:text-blue-700 md:p-0">Админ панель</a>
                          </li>
                        </Link>
                      :null}
                    </ul>
                  :null}
                </div>
                <div class="flex items-center lg:order-3">
                  <div class="flex-shrink-0 dropdown">
                    <a href="#" class="d-block link-body-emphasis text-decoration-none dropdown-toggle" data-bs-toggle="dropdown" aria-expanded="false">
                      <img src="https://github.com/mdo.png" alt="mdo" width="32" height="32" class="rounded-circle" />
                    </a>
                    <ul class="dropdown-menu text-small shadow">
                      {/*
                        {permissions?.map(per =>
                            <li><a class="dropdown-item" href="#">{per.topic} {per.all_records? `${per.all_records}`: 'безлимитно'}</a></li>
                        )}
                      */}
                    </ul>
                  </div>
                </div>
                <div class="flex items-center lg:order-2">
                  <button type="button" class="hidden sm:inline-flex items-center justify-center text-white bg-primary-700 hover:bg-primary-800 focus:ring-4 focus:ring-primary-300 font-medium rounded-lg text-xs px-3 py-1.5 mr-2 dark:bg-primary-600 dark:hover:bg-primary-700 focus:outline-none dark:focus:ring-primary-800"><svg aria-hidden="true" class="mr-1 -ml-1 w-5 h-5" fill="currentColor" viewBox="0 0 20 20" xmlns="http://www.w3.org/2000/svg"><path fill-rule="evenodd" d="M10 5a1 1 0 011 1v3h3a1 1 0 110 2h-3v3a1 1 0 11-2 0v-3H6a1 1 0 110-2h3V6a1 1 0 011-1z" clip-rule="evenodd"></path></svg> Действие</button>
                </div>
              </div>
            </nav>
            <nav class="bg-white border-gray-200 px-6">
              <div class="w-full">
                <div class="flex justify-between items-center">
                  <button id="toggleSidebar" aria-expanded="true" aria-controls="sidebar" class="hidden p-2 mr-3 text-gray-600 rounded cursor-pointer lg:inline hover:text-gray-900 hover:bg-gray-100 dark:text-gray-400 dark:hover:text-white dark:hover:bg-gray-700" data-bs-toggle="collapse" data-bs-target="#sidebar" aria-label="Toggle navigation">
                    <svg class="w-6 h-6" fill="currentColor" viewBox="0 0 20 20" xmlns="http://www.w3.org/2000/svg"><path fill-rule="evenodd" d="M3 5a1 1 0 011-1h12a1 1 0 110 2H4a1 1 0 01-1-1zM3 10a1 1 0 011-1h6a1 1 0 110 2H4a1 1 0 01-1-1zM3 15a1 1 0 011-1h12a1 1 0 110 2H4a1 1 0 01-1-1z" clip-rule="evenodd"></path></svg>
                  </button>
                </div>
              </div>
            </nav>
          </header >
          <main>
            <div>
              <div className="container-fluid">
                <div className="row align-items-stretch b-height">
                  <aside id="sidebar" className="h-screen col-md-2 my-3 bg-white collapse show width border rounded-3 g-0">
                    <div className="accordion accordion-flush" id="accordionFlushExample">
                        <div class="ml-5 grow items-center justify-between hidden w-full md:flex md:w-auto md:order-1" id="navbar-sticky">
                        </div>
                    </div>
                  </aside>
                  <section class="col p-3 m-3 border rounded-3 bg-white overflow-auto">
                    <div class="accordion accordion-flush" id="accordion">

                    </div>
                    <div>
                      <div class="bd-example">
                        <div class="tab-content" id="myTabContent">
                            <div class="container-fluid g-0">
                                <div class="relative mt-1 w-full">
                                  <div class="absolute inset-y-0 left-0 flex items-center pl-3 pointer-events-none">
                                    <svg class="w-4 h-4 text-gray-500 dark:text-gray-400" aria-hidden="true" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 20 20">
                                      <path stroke="currentColor" stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="m19 19-4-4m0-7A7 7 0 1 1 1 8a7 7 0 0 1 14 0Z" />
                                    </svg>
                                  </div>
                                  <input
                                    class="py-3 bg-gray-50 border border-gray-300 text-gray-900 sm:text-sm rounded-lg focus:ring-primary-500 focus:border-primary-500 block w-full pl-10 p-2.5"
                                    id="search"
                                    type="text"
                                    name="search_field"
                                    placeholder="Задайте вопрос"
                                    value={sendMessage}
                                    onChange={this.changeQueryText}
                                    aria-label="Search" />
                                  {/*<button type="submit" disabled={loading} value="Перевести" onClick={() => this.translateQuery()} class="text-white absolute right-2.5 bottom-2.5 bg-blue-700 hover:bg-blue-800 focus:ring-4 focus:outline-none focus:ring-blue-300 font-medium rounded-lg text-sm px-3 py-2">Перевести</button>*/}
                                  <button type="submit" disabled={loading} value="Отправить" onClick={() => this.getRequest()} class="text-white absolute right-2.5 bottom-2.5 bg-blue-700 hover:bg-blue-800 focus:ring-4 focus:outline-none focus:ring-blue-300 font-medium rounded-lg text-sm px-4 py-2">Отправить</button>
                                </div>
                            </div>
                            <div name='chat'>
                                {loading?
                                    <p>loading...</p>
                                :null}
                                <br />
                                {messages?.toReversed().map(m =>
                                    <>
                                        <div className="flex flex-row">
                                            <p>{m.request}</p>
                                        </div>
                                        <div className="flex flex-row-reverse">
                                            {m.status > 299?
                                                <p style={{ color: 'red' }}>{m.response}</p>
                                            :
                                                <p style={{ color: 'green' }}>{m.response}</p>
                                            }

                                        </div>
                                        <br />
                                    </>
                                )}
                            </div>
                        </div>
                      </div>
                    </div>
                  </section>
                </div>
              </div>
            </div>
          </main>
        </>
      )
    }
  }
}